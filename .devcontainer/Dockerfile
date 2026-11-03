FROM ghcr.io/acts-project/spack-container:24.5.0_linux-ubuntu26.04-x86_64_gcc-15.2.0_cxx23

RUN apt-get update && apt-get install -y clangd vim

RUN uv tool install -U --no-cache --no-build --force uv==0.11.24

RUN uv tool install --no-build --no-cache prek==0.4.5
