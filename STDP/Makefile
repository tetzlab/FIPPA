.DEFAULT_GOAL := comparison_brian2_arbor.png

BUILD_CATALOGUE_SCRIPT := build-catalogue

custom-catalogue.so: $(wildcard mechanisms/*.mod)
	$(BUILD_CATALOGUE_SCRIPT) custom mechanisms

.PRECIOUS: arbor_traces_%.dat arbor_spikes_%.dat
arbor_traces_%.dat arbor_spikes_%.dat: config_%.json arbor_lif_stdp.py custom-catalogue.so
	./run_arbor.sh $*

.PRECIOUS: brian2_traces_%.dat brian2_spikes_%.dat
brian2_traces_%.dat brian2_spikes_%.dat: config_%.json brian2_lif_stdp.py
	./run_brian2.sh $*

comparison_%.png: arbor_traces_%.dat arbor_spikes_%.dat brian2_traces_%.dat brian2_spikes_%.dat compare.py
	./compare.py $*

.PHONY: clean
clean:
	rm -f *png *dat *so

.PHONY: help
help:
	@echo "For configuration variant config_brian2_arbor.json do: make comparison_brian2_arbor.png"
