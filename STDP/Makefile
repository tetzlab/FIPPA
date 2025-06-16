.DEFAULT_GOAL := comparison_brian2_arbor.png

custom-catalogue.so: $(wildcard mechanisms/*.mod)
	./build_arbor_catalogue.sh custom mechanisms

.PRECIOUS: arbor_traces_%_lif.dat arbor_spikes_%_lif.dat arbor_traces_%_classical.dat
arbor_traces_%_lif.dat arbor_spikes_%_lif.dat: config_%_lif.json arbor_stdp_lif.py arbor_stdp_classical.py custom-catalogue.so config_brian2_arbor_lif.json config_brian2_arbor_classical.json
	./run_arbor.sh $*

.PRECIOUS: brian2_traces_%_lif.dat brian2_spikes_%_lif.dat brian2_traces_%_classical.dat
brian2_traces_%_lif.dat brian2_spikes_%_lif.dat brian2_traces_%_classical.dat: config_%_lif.json brian2_stdp_lif.py brian2_stdp_classical.py config_brian2_arbor_lif.json config_brian2_arbor_classical.json
	./run_brian2.sh $*

comparison_%.png: arbor_traces_%_lif.dat arbor_spikes_%_lif.dat brian2_traces_%_lif.dat brian2_spikes_%_lif.dat brian2_traces_%_classical.dat compare.py
	./compare.py $* > goodness_of_fit.txt

.PHONY: clean
clean:
	rm -f *png *dat *so

.PHONY: help
help:
	@echo "For configuration variant config_brian2_arbor.json do: make comparison_brian2_arbor.png"
