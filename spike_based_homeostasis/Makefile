.DEFAULT_GOAL := homeostasis.svg

BUILD_CATALOGUE_SCRIPT := arbor-build-catalogue
CONFIG := config_arbor.json

homeostasis-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) homeostasis mechanisms

arbor_spikes.dat arbor_traces.dat: homeostasis-catalogue.so $(CONFIG) arbor_homeostasis.py
	./arbor_homeostasis.py $(CONFIG) --catalogue ./$<

brian2_spikes.dat brian2_traces.dat: config_brian2.json brian2_homeostasis.py
	./brian2_homeostasis.py config_brian2.json

homeostasis.svg: arbor_input.dat arbor_traces.dat brian2_input.dat brian2_traces.dat compare.py
	./compare.py $(CONFIG) --save homeostasis.svg

.PHONY: clean
clean:
	rm -f *svg *so *dat
