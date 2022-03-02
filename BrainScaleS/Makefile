BUILD_CATALOGUE_SCRIPT := spack build-env arbor arbor-build-catalogue

.PHONY : all
all: plateau.png

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
	@echo "For configuration variant config_brian2_arbor.json do: make comparison_brian2_arbor.png"
