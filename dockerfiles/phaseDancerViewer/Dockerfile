FROM ubuntu:focal

RUN apt-get update -y
RUN apt-get upgrade -y 

RUN DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata

RUN apt-get install git -y

RUN apt-get install software-properties-common -y
RUN add-apt-repository ppa:deadsnakes/ppa -y
RUN apt-get update -y
RUN apt-get install python3.8 -y
RUN apt-get install python-is-python3 -y
RUN apt-get install python3-pip -y

################################ PHASEDANCER VIEWER ################################

RUN mkdir /phaseDancerViewer_bin
RUN mkdir /phaseDancerViewer_data
WORKDIR /phaseDancerViewer_bin
RUN git clone https://github.com/bposzewiecka/phaseDancerViewer.git
WORKDIR /phaseDancerViewer_bin/phaseDancerViewer
RUN pip install -r requirements.txt
RUN ln -s  /phaseDancerViewer_data/data static/data
RUN ln -s  /phaseDancerViewer_data/config.yaml static/config.yaml
RUN python manage.py migrate
RUN python manage.py shell < create_db_objects.py
ENV DJANGO_SUPERUSER_PASSWORD=johnny
RUN python manage.py createsuperuser --no-input --username admin --email johnny.bean@wp.pl

ENTRYPOINT ["python", "manage.py", "runserver", "0.0.0.0:8000"]
