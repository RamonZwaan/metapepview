FROM continuumio/miniconda3
RUN conda config --add channels bioconda && \
    conda install -y \
    python=3.11 \
    numpy=1.26 \
    pandas=2.1 \
    scikit-learn \
    scipy \
    pyteomics \
    matplotlib \
    seaborn \
    xlsxwriter && \
    conda install -c conda-forge dash=2.13.0 && \
    conda install -c conda-forge dash-bootstrap-components
COPY ./data /home/data
COPY ./src /home/src
WORKDIR /home
EXPOSE 8050
ENTRYPOINT python ./src/dashboard/index.py