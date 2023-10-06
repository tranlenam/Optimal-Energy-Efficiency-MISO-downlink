# hash:sha256:5b999ddd49d657f994ff2bd965f526cc98b26ce5a9160880d56b228d282452e7
FROM registry.codeocean.com/codeocean/matlab:2019a-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

COPY postInstall /
RUN /postInstall
