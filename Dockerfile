FROM continuumio/anaconda:5.2.0
MAINTAINER tmajaria@broadinstitute.org

RUN cd / && \
	git clone https://github.com/manning-lab/ldsc_wdl.git

RUN cd / && \
	git clone https://github.com/bulik/ldsc.git && \
	cd ldsc && \
	conda env create --file environment.yml

RUN [ "/bin/bash", "-c", "source activate ldsc" ]