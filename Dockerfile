FROM mambaorg/micromamba:1.5.6

USER root
RUN echo "user:x:1001:1001::/home/user:/bin/bash" >> /etc/passwd && \
    mkdir -p /home/user && chown -R 1001:1001 /home/user
RUN mkdir -p $HOME/.conda && chown -R 1000:1000 $HOME
COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes
RUN conda run --no-capture-output -n nf-ggi \
    python -c 'from rf2t_micro.weights import get_model_weights; get_model_weights(); from yunta.weights import get_model_weights; get_model_weights()'

ENV PATH=/opt/conda/envs/env/bin:$PATH
