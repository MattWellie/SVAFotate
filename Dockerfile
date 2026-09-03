# Minimal runtime image for svafotate_cpg. Everything installs from wheels, so no compiler is needed.
FROM python:3.13-slim

COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv

WORKDIR /app
COPY pyproject.toml README.md LICENSE uv.lock ./
COPY src ./src

RUN uv sync --no-dev --no-build --frozen --no-install-project && \
    uv pip install --no-deps . && \
    /app/.venv/bin/svafotate --version

ENV PATH="/app/.venv/bin:$PATH"
ENTRYPOINT ["svafotate"]
