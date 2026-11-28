"""QDP++ double-field HMC noise in global x-fast site order."""
import math
import numpy as np


MODULUS = 1 << 47
MULTIPLIER = 84000335758957
MASK = np.uint64(MODULUS - 1)
# U gamma_Chroma U^dagger = gamma_Bridge, in all four directions.
CHROMA_TO_BRIDGE = np.array([[0,1,0,1],[-1,0,-1,0],[0,-1,0,1],[1,0,-1,0]], dtype=float) / math.sqrt(2)

def seed_to_float(raw):
    # QDP++ rounds each multiply/add of the four base-4096 words to float.
    # A single float cast of the completed 47-bit fraction changes rare draws.
    value = (raw & np.uint64(4095)).astype(np.float32)
    for shift in (12, 24, 36):
        word = ((raw >> np.uint64(shift)) & np.uint64(4095)).astype(np.float32)
        value = value * np.float32(1 / 4096) + word
    return value * np.float32(1 / 2048)

def pack_seed(words):
    if len(words) != 4 or any(not isinstance(n, int) or not 0 <= n < limit
                              for n, limit in zip(words, (4096, 4096, 4096, 2048))):
        raise ValueError("invalid four-word QDP seed")
    seed = sum(word << (12 * i) for i, word in enumerate(words))
    if not seed:
        raise ValueError("a zero QDP seed cannot generate Gaussian noise")
    return seed


class QDPNoise:
    """Supply the native generator's real/imaginary draw pairs and conditional uniform."""
    def __init__(self, dims, words):
        self.seed = pack_seed(words)
        self.volume = math.prod(dims)
        self.stride = pow(MULTIPLIER, self.volume + 1, MODULUS)
        # The lower 47 bits of a uint64 product give the product modulo 2^47.
        powers = np.full(self.volume, MULTIPLIER, dtype=np.uint64)
        self.site_powers = np.bitwise_and(np.multiply.accumulate(powers), MASK)
        sites = np.arange(self.volume)
        parity = np.zeros(self.volume, dtype=np.uint8)
        for extent in dims:
            parity ^= (sites % extent % 2).astype(np.uint8)
            sites //= extent
        self.odd_powers = self.site_powers[parity.astype(bool)]
        self.imaginary = None

    @property
    def words(self):
        if self.imaginary is not None:
            raise ValueError("incomplete QDP Gaussian draw")
        return [(self.seed >> (12*i)) & 4095 for i in range(4)]

    def _uniforms(self, powers, seed):
        raw = (powers * np.uint64(seed)) & MASK
        # sranf returns float even when the destination field is double.
        return seed_to_float(raw).astype(np.float64)

    def _gaussian(self, powers, components):
        result = np.empty((len(powers), components), dtype=np.complex128)
        jump = pow(self.stride, 2 * components, MODULUS)
        radius_seed = self.seed
        angle_seed = self.seed * jump % MODULUS
        two_fields = self.stride * self.stride % MODULUS
        for component in range(components):
            radius = np.sqrt(-2.0 * np.log(self._uniforms(powers, radius_seed)))
            # This decimal constant is used by QDP++ fill_gaussian(RComplex).
            angle = 6.283185307 * self._uniforms(powers, angle_seed)
            result[:, component].real = radius * np.cos(angle)
            result[:, component].imag = radius * np.sin(angle)
            radius_seed = radius_seed * two_fields % MODULUS
            angle_seed = angle_seed * two_fields % MODULUS
        self.seed = angle_seed
        return result

    def standard_normal(self, shape):
        shape = tuple(shape)
        if self.imaginary is not None:
            if self.imaginary.shape != shape:
                raise ValueError("inconsistent real and imaginary QDP draw shapes")
            part, self.imaginary = self.imaginary, None
            return part
        if shape == (self.volume, 4, 3, 3):
            field = np.empty(shape, dtype=np.complex128)
            for mu in range(4):
                field[:, mu] = self._gaussian(self.site_powers, 9).reshape(-1,3,3)
        elif shape == (len(self.odd_powers), 4, 3):
            chroma = self._gaussian(self.odd_powers, 12).reshape(shape)
            field = np.einsum("ab,vbc->vac", CHROMA_TO_BRIDGE, chroma)
        else:
            raise ValueError(f"unsupported QDP HMC draw shape {shape}")
        self.imaginary = field.imag
        return field.real

    def random(self):
        if self.imaginary is not None:
            raise ValueError("Metropolis draw before completing a Gaussian field")
        self.seed = self.seed * MULTIPLIER % MODULUS
        return float(seed_to_float(np.uint64(self.seed)))
