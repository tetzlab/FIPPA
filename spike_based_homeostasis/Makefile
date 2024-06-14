.DEFAULT_GOAL := homeostasis.svg

BUILD_CATALOGUE_SCRIPT := arbor-build-catalogue
CONFIG := config_arbor.json
NUM_TRIALS := 100

homeostasis-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) homeostasis mechanisms

arbor_input.dat arbor_traces.dat arbor_spikes.dat: homeostasis-catalogue.so $(CONFIG) arbor_homeostasis.py
	./run_arbor.sh $(CONFIG) $(NUM_TRIALS) --catalogue ./$<

brian2_input.dat brian2_traces.dat brian2_spikes.dat: config_brian2.json brian2_homeostasis.py
	./run_brian2.sh config_brian2.json $(NUM_TRIALS)

homeostasis.svg: arbor_input.dat arbor_traces.dat arbor_spikes.dat brian2_input.dat brian2_traces.dat brian2_spikes.dat compare.py
	./compare.py $(CONFIG) $(NUM_TRIALS) --save homeostasis.svg

.PHONY: clean
clean:
	rm -f *svg *so *dat
