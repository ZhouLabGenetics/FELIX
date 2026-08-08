FROM ubuntu:20.04

LABEL maintainer="lhu@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get -y install --no-install-recommends \
      curl \
      ca-certificates \
      g++ \
      libhts-dev \
      make \
      pkg-config \
      python3 \
      tabix \
      time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'

RUN curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx && \
    mv plink-ng-2.0.0-a.6.16 plink-ng && \
    pixi run bash -c 'CC_BIN=$(ls .pixi/envs/default/bin/*-conda-linux-gnu-cc | head -1); "$CC_BIN" -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o plink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate' && \
    mv plink2_includes.a .pixi/envs/default/lib

RUN pixi run R CMD INSTALL .

# Build the FELIXla command (single binary; the packer, the RFMix/TRACTOR
# converters, export and query all dispatch from it) and run its smoke test.
#
# It is linked against the pixi environment's htslib, not the distro one:
# focal ships htslib 1.10, which predates hts_idx_nseq(). The rpath is the
# fixed in-image env path, so the binary resolves libhts.so.3 at runtime.
# CXXFLAGS is passed explicitly because the Makefile's default includes
# -march=native, which would tie the image to the build machine's CPU.
RUN mkdir -p /scripts && \
   cp -R tools/felixla/Makefile tools/felixla/src tools/felixla/tests tools/felixla/testdata tools/felixla/resources /scripts/ && \
   pixi run make -C /scripts CXX=g++ \
       CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic" \
       HTSLIB_CFLAGS="-I/app/.pixi/envs/default/include" \
       HTSLIB_LIBS="-L/app/.pixi/envs/default/lib -lhts -Wl,-rpath,/app/.pixi/envs/default/lib" && \
   pixi run make -C /scripts test CXX=g++ \
       CXXFLAGS="-O3 -std=c++17 -Wall -Wextra -Wpedantic" \
       HTSLIB_CFLAGS="-I/app/.pixi/envs/default/include" \
       HTSLIB_LIBS="-L/app/.pixi/envs/default/lib -lhts -Wl,-rpath,/app/.pixi/envs/default/lib" && \
   /scripts/bin/felixla --help | grep -q "felixla \[input flags\]"

ENV FELIXLA_HOME=/scripts
ENV FELIXLA_RESOURCES=/scripts/resources

# Expose the complete pixi environment (including Python + NumPy) so the
# tutorial runner can generate its simulated data inside Docker/Singularity.
ENV PATH=/scripts/bin:/app/.pixi/envs/default/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/ && \
   chmod a+x /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step2_SPAtests.R /usr/local/bin/step3_LDmat.R /usr/local/bin/createSparseGRM.R && \
   createSparseGRM.R --help && \
   step1_fitNULLGLMM.R --help && \
   step2_SPAtests.R --help && \
   step3_LDmat.R --help
