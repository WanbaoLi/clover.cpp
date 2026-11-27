#!/bin/sh

echo ""
echo "################################################################################"
echo "######################## Bridge++ Makefile Setup script ########################"
echo "################################################################################"
#
echo ""

datestr=`date +%Y%m%d%H%M`
if [ -f Makefile ]; then
    echo "Current Makefile is saved to Makefile.$datestr."
    cp -p Makefile Makefile.$datestr
    echo ""
fi

#----------------------------------------------------------------
echo "### Select Platform ###"
echo "1 [PC WorkStation with GNU C++] (default)"
echo "2 [PC WorkStation with Clang C++]"
echo "3 [PC WorkStation with Intel C++]"
echo "4 [PC WorkStation with NVIDIA HPC SDK]"
echo "5 [Fugaku_CLANG]"
echo "6 [SX_Aurora@KEK]"
echo ""

echo -n "Platform? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	TARGET_FLAG="PC_GNU"
	echo "You have selected [PC Workstation with GNU C++]"
	;;
    2)
	TARGET_FLAG="PC_CLANG"
	echo "You have selected [PC Workstation with Clang C++]"
	;;
    3)
	TARGET_FLAG="PC_INTEL"
	echo "You have selected [PC Workstation with Intel C++]"
	;;
    4)
	TARGET_FLAG="PC_NVIDIA"
	echo "You have selected [PC Workstation with NVIDIA HPC SDK]"
	;;
    5)
	TARGET_FLAG="Fugaku_CLANG"
	echo "You have selected [Fugaku_CLANG]"
	;;
    6)
	TARGET_FLAG="SX_AURORA"
	echo "You have selected [SX_AURORA]"
	;;
    *)
	TARGET_FLAG="PC_GNU"
	echo "You have selected default value [PC Workstation With GNU C++]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Communicator Type ###"
echo "1 [Single] (default)"
echo "2 [MPI]"
echo "3 [FJMPI (for Fugkau, MPI with Fujitus extensnion)]"
echo ""

echo -n "Communication Type? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_MPI_FLAG="no"
	echo "You have selected [Single]"
	;;
    2)
	USE_MPI_FLAG="yes"
	echo "You have selected [MPI]"
	;;
    3)
	USE_MPI_FLAG="fjmpi"
	echo "You have selected [FJMPI]"
	;;
    *)
	USE_MPI_FLAG="no"
	echo "You have selected default value [Single]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Multi-threading Option ###"
echo "1 [OpenMP]"
echo "2 [No] (default)"
echo ""

echo -n "Use Multi-threading? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_THREAD_FLAG="omp"
	echo "You have selected [OpenMP]"
	;;
    2)
	USE_THREAD_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_THREAD_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Optimized Code ###"

if [ "${USE_THREAD_FLAG}" != "no" ]; then

    echo "## Optimized code is required when multi-threading is enabled."
    echo
    echo "[Yes] enforced."

    USE_OPT_CODE_FLAG="yes"

else
    echo "1 [Yes] (default)"
    echo "2 [No]"
    echo ""

    echo -n "Use Optimized Code? > "
    read VALUE

    echo ""

    case ${VALUE} in
	1)
	    USE_OPT_CODE_FLAG="yes"
	    echo "You have selected [Yes]"
	    ;;
	2)
	    USE_OPT_CODE_FLAG="no"
	    echo "You have selected [No]"
	    ;;
	*)
	    USE_OPT_CODE_FLAG="yes"
	    echo "You have selected default value [Yes]"
	    ;;
    esac
fi

echo ""

#----------------------------------------------------------------
echo "### Select Compiler Optimization Level ###"
echo "1 [Low]"
echo "2 [Standard] (default)"
echo "3 [High]"
echo "4 [Trace] only for SX-Aurora"
echo ""

echo -n "Optimization Level? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_OPT_LEVEL_FLAG="low"
	echo "You have selected [Low]"
	;;
    2)
	USE_OPT_LEVEL_FLAG="std"
	echo "You have selected [Standard]"
	;;
    3)
	USE_OPT_LEVEL_FLAG="high"
	echo "You have selected [High]"
	;;
    4)
	USE_OPT_LEVEL_FLAG="trace"
	echo "You have selected [Trace]"
	;;
    *)
	USE_OPT_LEVEL_FLAG="std"
	echo "You have selected default value [Standard]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Static linkage Option ###"
echo "## [Yes] needs static system libraries"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Static linkage? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_STATIC_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_STATIC_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_STATIC_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Debug Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Debug Option? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_DEBUG_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_DEBUG_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_DEBUG_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Select Complex Type ###"
echo "1 [STL complex] (default)"
echo "2 [C99 complex]"
echo ""

echo -n "Complex Type? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_COMPLEX_FLAG="std"
	echo "You have selected [STL complex]"
	;;
    2)
	USE_COMPLEX_FLAG="c99"
	echo "You have selected [C99 complex]"
	;;
    *)
	USE_COMPLEX_FLAG="std"
	echo "You have selected default value [STL complex]"
	;;
esac

echo ""

echo "### Select Use C++11 Functions ###"
echo "1 [Yes] (default)"
echo "2 [No]"
echo ""

echo -n "Use cpp11? > "
read VALUE

echo ""

case ${VALUE} in
    1)
        USE_CPP11_FLAG="yes"
        echo "You have selected [Yes]"
        ;;
    2)
        USE_CPP11_FLAG="no"
        ;;
    *)
        USE_CPP11_FLAG="yes"
        echo "You have selected default value [yes]"
        ;;
esac

echo ""

#----------------------------------------------------------------
echo "### Select Enable Python Interface ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use python? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_PYTHON_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_PYTHON_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_PYTHON_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo

#----------------------------------------------------------------
echo "### Use Alternative Code ###"
echo "##   architecture specific alternative implementations"
echo "##   set yes for QXS (Fugaku) or Accel (GPU)"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Alternative Code? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_ALT_CODE_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_ALT_CODE_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_ALT_CODE_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo

#----------------------------------------------------------------
echo "### Use Alternative QXS library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Alternative QXS library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_ALT_QXS_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_ALT_QXS_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_ALT_QXS_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
if [ "${USE_ALT_QXS_FLAG}" = "yes" ]; then

  echo "### Use QXS Architecture Option ###"
  echo "## acle is for Fugaku"
  echo "1 [general] (default)"
  echo "2 [acle]"
  echo ""

  echo -n "Use QXS Architecture? > "
  read VALUE

  echo ""

  case ${VALUE} in
    1)
      USE_QXS_ARCH_FLAG="general"
      echo "You have selected [general]"
      ;;
    2)
      USE_QXS_ARCH_FLAG="acle"
      echo "You have selected [acle]"
      ;;
    *)
      USE_QXS_ARCH_FLAG="general"
      echo "You have selected default value [general]"
      ;;
  esac
else
  USE_QXS_ARCH_FLAG="general"
fi
echo ""


#----------------------------------------------------------------
echo "### Use Alternative Accel library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use Alternative Accel library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_ALT_ACCEL_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_ALT_ACCEL_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_ALT_ACCEL_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
if [ "${USE_ALT_ACCEL_FLAG}" = "yes" ]; then

  echo "### Use Accel Implementation Option ###"
  echo "## currently only openacc is avaiable"
  echo "1 [openacc] (default)"
  echo ""

  echo -n "Use Accel Implementation? > "
  read VALUE

  echo ""

  case ${VALUE} in
    1)
      USE_ACCEL_IMPL_FLAG="openacc"
      echo "You have selected [openacc]"
      ;;
    *)
      USE_ACCEL_IMPL_FLAG="openacc"
      echo "You have selected default value [openacc]"
      ;;
  esac
  echo

#-------------

  echo "### Use Device Architecture Option ###"
  echo "## currently only openacc is avaiable"
  echo "## auto: using nvidia-smi, which might not be available in some environments"
  echo "##       (such as on systems at Univ. Tsukuba)."
  echo "1 [kepler]"
  echo "2 [pascal]"
  echo "3 [volta]"
  echo "4 [ampere]"
  echo "5 [hopper]"
  echo "6 [rumors]"
  echo "7 [auto] (default)"
  echo ""

  echo -n "Use Device Architecture Option? > "
  read VALUE

  echo ""

  case ${VALUE} in
    1)
      USE_DEVICE_ARCH_FLAG="kepler"
      echo "You have selected [kepler]"
      ;;
    2)
      USE_DEVICE_ARCH_FLAG="pascal"
      echo "You have selected [pascal]"
      ;;
    3)
      USE_DEVICE_ARCH_FLAG="volta"
      echo "You have selected [volta]"
      ;;
    4)
      USE_DEVICE_ARCH_FLAG="ampere"
      echo "You have selected [ampere]"
      ;;
    5)
      USE_DEVICE_ARCH_FLAG="hopper"
      echo "You have selected [hopper]"
      ;;
    6)
      USE_DEVICE_ARCH_FLAG="rumors"
      echo "You have selected [rumors]"
      ;;
    7)
      USE_DEVICE_ARCH_FLAG="auto"
      echo "You have selected [auto]"
      ;;
    *)
      USE_DEVICE_ARCH_FLAG="auto"
      echo "You have selected default value [auto]"
      ;;
  esac
else
  USE_ACCEL_IMPL_FLAG="openacc"
  USE_DEVICE_ARCH_FLAG="auto"
fi

echo ""



#----------------------------------------------------------------
echo "### Build directory ###"
echo
echo "which directory to build and store library and include files?"
echo -n "[default: build] > "
read VALUE

if [ -n "${VALUE}" ]; then
    BUILD_DIR=${VALUE}
else
    BUILD_DIR="build"
fi

echo

#----------------------------------------------------------------
echo "### Build Tests directory ###"
echo
echo "which directory to build Tests and store executables?"
echo -n "[default: ${BUILD_DIR}/tests] > "
read VALUE

if [ -n "${VALUE}" ]; then
    TEST_DIR=${VALUE}
else
    TEST_DIR="${BUILD_DIR}/tests"
fi

echo

#----------------------------------------------------------------
echo "### Select Gauge group ###"
echo "1 [SU(3)] (default)"
echo "2 [SU(2)]"
echo "3 [General SU(N)]"
echo ""

echo -n "Gauge group? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_GAUGE_GROUP_FLAG="su3"
	echo "You have selected [SU(3)]"
	;;
    2)
	USE_GAUGE_GROUP_FLAG="su2"
	echo "You have selected [SU(2)]"
	;;
    3)
	USE_GAUGE_GROUP_FLAG="general"
	echo "You have selected [General SU(N)]"
	;;
    *)
	USE_GAUGE_GROUP_FLAG="su3"
	echo "You have selected default value [SU(3)]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Select Target test name ###"
echo "1 [all] (default)"
echo "2 [Eigensolver]"
echo "3 [FFT]"
echo "4 [Gauge]"
echo "5 [GradientFlow]"
echo "6 [HMC]"
echo "7 [HotStart]"
echo "8 [IO]"
echo "9 [Mult]"
echo "10 [PolyakovLoop]"
echo "11 [QuarkNumberSusceptibility]"
echo "12 [RandomNumbers]"
echo "13 [Rational]"
echo "14 [SF_fAfP]"
echo "15 [Solver]"
echo "16 [Spectrum]"
echo "17 [TopologicalCharge]"
echo "18 [WilsonLoop]"
echo "19 [none]"
echo ""

echo -n "Test name? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	TEST_FLAG="all"
	echo "You have selected [all]"
	;;
    2)
	TEST_FLAG="Eigensolver"
	echo "You have selected [Eigensolver]"
	;;
    3)
	TEST_FLAG="FFT"
	echo "You have selected [FFT]"
	;;
    4)
	TEST_FLAG="Gauge"
	echo "You have selected [Gauge]"
	;;
    5)
	TEST_FLAG="GradientFlow"
	echo "You have selected [GradientFlow]"
	;;
    6)
	TEST_FLAG="HMC"
	echo "You have selected [HMC]"
	;;
    7)
	TEST_FLAG="HotStart"
	echo "You have selected [HotStart]"
	;;
    8)
	TEST_FLAG="IO"
	echo "You have selected [IO]"
	;;
    9)
	TEST_FLAG="Mult"
	echo "You have selected [Mult]"
	;;
    10)
	TEST_FLAG="PolyakovLoop"
	echo "You have selected [PolyakovLoop]"
	;;
    11)
	TEST_FLAG="QuarkNumberSusceptibility"
	echo "You have selected [QuarkNumberSusceptibility]"
	;;
    12)
	TEST_FLAG="RandomNumbers"
	echo "You have selected [RandomNumbers]"
	;;
    13)
	TEST_FLAG="Rational"
	echo "You have selected [Rational]"
	;;
    14)
	TEST_FLAG="SF_fAfP"
	echo "You have selected [SF_fAfP]"
	;;
    15)
	TEST_FLAG="Solver"
	echo "You have selected [Solver]"
	;;
    16)
	TEST_FLAG="Spectrum"
	echo "You have selected [Spectrum]"
	;;
    17)
	TEST_FLAG="TopologicalCharge"
	echo "You have selected [TopologicalCharge]"
	;;
    18)
	TEST_FLAG="WilsonLoop"
	echo "You have selected [WilsonLoop]"
	;;
    19)
	TEST_FLAG="none"
	echo "You have selected [none]"
	;;
    *)
	TEST_FLAG="all"
	echo "You have selected default value [all]"
	;;
esac

echo ""
#----------------------------------------------------------------
echo "### Use Test Manager ###"
echo "1 [Yes] (default)"
echo "2 [No]"
echo ""

echo -n "Use Test Manager? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_TESTMANAGER_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_TESTMANAGER_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_TESTMANAGER_FLAG="yes"
	echo "You have selected default value [Yes]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### External Library Options ###"
echo ""

#----------------------------------------------------------------
echo "### Extra library directory ###"
echo
echo "which directory to locate extra libraries?"
echo -n "[default: extra] > "
read VALUE

if [ -n "${VALUE}" ]; then
    EXTRA_DIR=${VALUE}
else
    EXTRA_DIR="extra"
fi

echo

#----------------------------------------------------------------
echo "### Use LIME library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use LIME library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_LIME_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_LIME_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_LIME_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use SFMT(SIMD Fast Mersenne Twister) library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use SFMT library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_SFMT_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_SFMT_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_SFMT_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use NTL library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use NTL library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_NTL_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_NTL_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_NTL_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use FFTW library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use FFTW library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_FFTW_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_FFTW_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_FFTW_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use YAML-CPP library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use YAML-CPP library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_YAMLCPP_LINBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_YAMLCPP_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_YAMLCPP_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use tinyxml2 library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use tinyxml2 library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_TINYXML2_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_TINYXML2_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_TINYXML2_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""

#----------------------------------------------------------------
echo "### Use QWS library Option ###"
echo "1 [Yes]"
echo "2 [No] (default)"
echo ""

echo -n "Use QWS library? > "
read VALUE

echo ""

case ${VALUE} in
    1)
	USE_QWS_LIBRARY_FLAG="yes"
	echo "You have selected [Yes]"
	;;
    2)
	USE_QWS_LIBRARY_FLAG="no"
	echo "You have selected [No]"
	;;
    *)
	USE_QWS_LIBRARY_FLAG="no"
	echo "You have selected default value [No]"
	;;
esac

echo ""


#----------------------------------------------------------------



sed \
    -e "s/@@TARGET@@/$TARGET_FLAG/" \
    -e "s/@@USE_MPI@@/$USE_MPI_FLAG/" \
    -e "s/@@USE_THREAD@@/$USE_THREAD_FLAG/" \
    -e "s/@@USE_OPT_CODE@@/$USE_OPT_CODE_FLAG/" \
    -e "s/@@USE_OPT_LEVEL@@/$USE_OPT_LEVEL_FLAG/" \
    -e "s/@@USE_STATIC@@/$USE_STATIC_FLAG/" \
    -e "s/@@USE_DEBUG@@/$USE_DEBUG_FLAG/" \
    -e "s/@@USE_COMPLEX@@/$USE_COMPLEX_FLAG/" \
    -e "s/@@USE_CPP11@@/$USE_CPP11_FLAG/" \
    -e "s/@@USE_PYTHON@@/$USE_PYTHON_FLAG/" \
    -e "s/@@USE_ALT_CODE@@/$USE_ALT_CODE_FLAG/" \
    -e "s/@@USE_ALT_QXS@@/$USE_ALT_QXS_FLAG/" \
    -e "s/@@USE_QXS_ARCH@@/$USE_QXS_ARCH_FLAG/" \
    -e "s/@@USE_ALT_ACCEL@@/$USE_ALT_ACCEL_FLAG/" \
    -e "s/@@USE_ACCEL_IMPL@@/$USE_ACCEL_IMPL_FLAG/" \
    -e "s/@@USE_DEVICE_ARCH@@/$USE_DEVICE_ARCH_FLAG/" \
    -e "s/@@USE_GAUGE_GROUP@@/$USE_GAUGE_GROUP_FLAG/" \
    -e "s/@@TEST@@/$TEST_FLAG/" \
    -e "s/@@USE_TESTMANAGER@@/$USE_TESTMANAGER_FLAG/" \
    -e "s/@@USE_LIME_LIBRARY@@/$USE_LIME_LIBRARY_FLAG/" \
    -e "s/@@USE_SFMT_LIBRARY@@/$USE_SFMT_LIBRARY_FLAG/" \
    -e "s/@@USE_NTL_LIBRARY@@/$USE_NTL_LIBRARY_FLAG/" \
    -e "s/@@USE_FFTW_LIBRARY@@/$USE_FFTW_LIBRARY_FLAG/" \
    -e "s/@@USE_YAMLCPP_LIBRARY@@/$USE_YAMLCPP_LIBRARY_FLAG/" \
    -e "s/@@USE_TINYXML2_LIBRARY@@/$USE_TINYXML2_LIBRARY_FLAG/" \
    -e "s/@@USE_QWS_LIBRARY@@/$USE_QWS_LIBRARY_FLAG/" \
    -e "s,@@BUILD_DIR@@,$BUILD_DIR," \
    -e "s,@@TEST_DIR@@,$TEST_DIR," \
    -e "s,@@EXTRA_DIR@@,$EXTRA_DIR," \
Makefile.in > Makefile

#----------------------------------------------------------------
echo "Makefile is generated."

exit 0

