# Build context is staged by benchmark-run; never build against the repository root.
# This platform-specific manifest is also recorded in runtime.json.
FROM --platform=linux/amd64 python:3.11.11-slim-bookworm@sha256:a8e0a3090316aed0b11037aac613aef32fb1747dcc1dcb5c0f6c727a0113a07f AS hyphy-build
RUN rm -f /etc/apt/sources.list.d/debian.sources \
    && printf 'deb [check-valid-until=no] https://snapshot.debian.org/archive/debian/20250501T000000Z bookworm main\n' > /etc/apt/sources.list \
    && apt-get update && apt-get install -y --no-install-recommends g++ cmake make \
    && rm -rf /var/lib/apt/lists/*
ADD --checksum=sha256:59c4a7420f6a6f73d78f55ca13fc368b0893582fde468f5a48d25628b28177dc https://codeload.github.com/veg/hyphy/tar.gz/c2daaafe3f372e8e0f44e275db61f37c74a1516d /tmp/hyphy.tar.gz
RUN mkdir /tmp/hyphy && tar xzf /tmp/hyphy.tar.gz -C /tmp/hyphy --strip-components=1 \
    && cmake -S /tmp/hyphy -B /tmp/hyphy/build -DCMAKE_INSTALL_PREFIX=/opt/hyphy \
       -DCMAKE_DISABLE_FIND_PACKAGE_OpenMP=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_CURL=TRUE \
       -DCMAKE_DISABLE_FIND_PACKAGE_MPI=TRUE -DNOZLIB=ON -DNOBLAS=ON -DNOAVX=ON \
       -DNOSSE4=ON -DNONEON=ON -DTARGET_NATIVE=OFF \
       -DCMAKE_EXE_LINKER_FLAGS="-static-libgcc -static-libstdc++" \
    && cmake --build /tmp/hyphy/build --target hyphy --parallel 4 \
    && cmake --install /tmp/hyphy/build \
    && dpkg-query -W > /opt/hyphy/build-packages.txt

FROM --platform=linux/amd64 python:3.11.11-slim-bookworm@sha256:a8e0a3090316aed0b11037aac613aef32fb1747dcc1dcb5c0f6c727a0113a07f

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONHASHSEED=0 \
    PYTHONPATH=/app/src \
    OMP_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    HF_HUB_OFFLINE=1 \
    TRANSFORMERS_OFFLINE=1 \
    HF_DATASETS_OFFLINE=1 \
    HF_HOME=/tmp/huggingface \
    MPLCONFIGDIR=/tmp/matplotlib \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1

COPY --from=hyphy-build /opt/hyphy /opt/hyphy
ENV PATH=/opt/hyphy/bin:$PATH
RUN hyphy --version

WORKDIR /app
COPY requirements.lock /app/requirements.lock
RUN python -m pip install --index-url https://pypi.org/simple \
        --only-binary=:all: --require-hashes -r /app/requirements.lock \
    && python -m pip check \
    && python -c "import importlib.metadata, shutil; assert importlib.metadata.version('tn93') == '1.2.2'; assert shutil.which('tn93') is None, 'Unpinned external TN93 binary'"

# src contains only allowlisted production Python modules. No docs, data,
# manuscript, empirical reference, tests, agents or local model cache.
COPY src/ /app/src/
RUN find /app/src -type f ! -name '*.py' -print -quit | \
        python -c "import sys; assert not sys.stdin.read(), 'Non-Python source in inference snapshot'"

# Docker invocation supplies the host UID/GID, network=none, a read-only root
# filesystem, a /tmp tmpfs and only /inputs (read-only) plus /outputs mounts.
# No ENTRYPOINT: the sealed protocol records the explicit worker command.
