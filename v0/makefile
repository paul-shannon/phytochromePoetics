ID = 16may2020
IMAGE_NAME = phytochrome_poetics_$(ID)

default:
	@echo targets:  shiny build

build:
	docker build -t $(IMAGE_NAME) -f Dockerfile .


bash:
	docker run -it \
           -p 3838:9012 \
           --entrypoint /bin/bash \
           $(IMAGE_NAME)
