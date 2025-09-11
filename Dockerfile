FROM python:3

WORKDIR /home

COPY ./src ./src
COPY ./pyproject.toml ./pyproject.toml

ENV HOME="/home"
RUN pip install .

EXPOSE 8050
CMD [ "metapepview" ]