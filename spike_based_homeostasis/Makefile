.DEFAULT_GOAL := homeostasis.svg

BUILD_CATALOGUE_SCRIPT := arbor-build-catalogue
CONFIG := config.json

homeostasis-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) homeostasis mechanisms

spikes.dat traces.dat: homeostasis-catalogue.so $(CONFIG)
	./arbor_homeostasis.py $(CONFIG) --catalogue ./$<

homeostasis.svg: spikes.dat traces.dat
	./plot.py $(CONFIG) --save homeostasis.svg

.PHONY: clean
clean:
	rm -f *svg *so *dat
