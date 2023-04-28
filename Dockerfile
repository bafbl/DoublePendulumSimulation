FROM python:3.8-buster

RUN apt-get update && apt-get upgrade -y
#RUN apt-get install -y python3-pandas

RUN python3 -m pip --no-input install wheel openpyxl pandas


COPY main.py /

ENTRYPOINT ["/main.py"]
