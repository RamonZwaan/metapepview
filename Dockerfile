FROM python:3

WORKDIR /home

COPY ./requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY ./data ./data
COPY ./src ./src

EXPOSE 8050
CMD [ "python", "./src/index.py" ]