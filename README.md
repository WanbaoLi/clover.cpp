# Wilson–Clover CPU Base

本提交固定 Bridge++ 2.1.3 的裁剪库，提供 MPI、OpenMP、C-LIME、通用 CPU 和 512 位 ARM SVE2 构建。官方归档保存在父提交及 main 分支。此基线的 src 源自历史裁剪快照 86b0fc2，编译所需的 C-LIME 已源码化。

```sh
cmake -S . -B build -DCLOVER_QXS=general
cmake --build build --parallel 8
ctest --test-dir build --output-on-failure
cmake --install build --prefix "$PWD/install"
```

在 512 位 SVE2 计算节点上配置 `-DCLOVER_QXS=sve512`。自测使用时间反周期平面波，比较 Wilson–Clover 作用后的范数与独立解析式。编译需要 Linux、MPI、OpenMP 和支持 C++11 的工具链。
