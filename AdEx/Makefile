BUILD_CATALOGUE_SCRIPT := arbor-build-catalogue

.PHONY : all
all: tonic_spiking.png adaptation.png initial_burst.png regular_bursting.png delayed_accelerating.png delayed_regular_bursting.png transient_spiking.png irregular_spiking.png

custom-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) custom mechanisms

.PRECIOUS: arbor_traces_%.dat arbor_spikes_%.dat
arbor_traces_%.dat: config_%.json arbor_adex.py custom-catalogue.so
	./arbor_adex.py $*

%.png: arbor_traces_%.dat plot.py
	./plot.py $*

.PHONY: clean
clean:
	rm -f *png *dat *so

.PHONY: help
help:
	@echo "For configuration variant config_tonic_spiking.json do: make tonic_spiking.png"
