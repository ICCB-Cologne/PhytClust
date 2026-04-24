FROM python:3.11-slim

WORKDIR /app

# System deps for biopython/numpy/scipy
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Install Python deps first for better layer caching
COPY pyproject.toml .
RUN pip install --no-cache-dir \
    fastapi \
    "uvicorn[standard]" \
    jinja2 \
    python-multipart \
    numpy \
    biopython \
    matplotlib \
    pandas \
    scipy \
    PyYAML \
    rich

# Install package
COPY src/ src/
RUN pip install --no-cache-dir -e .

EXPOSE 8080

ENV PHYTCLUST_PUBLIC_MODE=1

CMD ["uvicorn", "phytclust.gui.api:app", "--host", "0.0.0.0", "--port", "8080"]
