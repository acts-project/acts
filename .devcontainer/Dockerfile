FROM ghcr.io/acts-project/spack-container:23.2.0_linux-ubuntu24.04_gcc-13.3.0_cxx20

RUN apt-get update && apt-get install -y clangd vim

RUN uv tool install --no-cache prek
