.DEFAULT_GOAL := comparison_brian2_arbor.png

custom-catalogue.so: $(wildcard mechanisms/*.mod)
	./build_arbor_catalogue.sh custom mechanisms

.PRECIOUS: arbor_traces_brian2_arbor_lif.dat arbor_spikes_brian2_arbor_lif.dat arbor_traces_brian2_arbor_classical.dat
arbor_traces_brian2_arbor_lif.dat arbor_spikes_brian2_arbor_lif.dat: arbor_stdp_lif.py arbor_stdp_classical.py custom-catalogue.so config_lif.json config_classical.json
	./run_arbor.sh "brian2_arbor"

.PRECIOUS: arbor_traces_nest_arbor_lif.dat arbor_spikes_nest_arbor_lif.dat arbor_traces_nest_arbor_classical.dat
arbor_traces_nest_arbor_lif.dat arbor_spikes_nest_arbor_lif.dat: arbor_stdp_lif.py arbor_stdp_classical.py custom-catalogue.so config_lif.json config_classical.json
	./run_arbor.sh "nest_arbor"

.PRECIOUS: brian2_traces_brian2_arbor_lif.dat brian2_spikes_brian2_arbor_lif.dat brian2_traces_brian2_arbor_classical.dat
brian2_traces_brian2_arbor_lif.dat brian2_spikes_brian2_arbor_lif.dat brian2_traces_brian2_arbor_classical.dat: brian2_stdp_lif.py brian2_stdp_classical.py config_lif.json config_classical.json
	./run_brian2.sh "brian2_arbor"

.PRECIOUS: nest_traces_nest_arbor_lif.dat nest_spikes_nest_arbor_lif.dat nest_traces_nest_arbor_classical.dat
nest_traces_nest_arbor_lif.dat nest_spikes_nest_arbor_lif.dat nest_traces_nest_arbor_classical.dat: nest_stdp_lif.py nest_stdp_classical.py config_lif.json config_classical.json
	./run_nest.sh "nest_arbor"

comparison_brian2_arbor.png: arbor_traces_brian2_arbor_lif.dat arbor_spikes_brian2_arbor_lif.dat brian2_traces_brian2_arbor_lif.dat brian2_spikes_brian2_arbor_lif.dat brian2_traces_brian2_arbor_classical.dat compare.py
	./compare.py "brian2_arbor" > goodness_of_fit.txt

comparison_nest_arbor.png: arbor_traces_nest_arbor_lif.dat arbor_spikes_nest_arbor_lif.dat nest_traces_nest_arbor_lif.dat nest_spikes_nest_arbor_lif.dat nest_traces_nest_arbor_classical.dat compare.py
	./compare.py "nest_arbor" > goodness_of_fit.txt

.PHONY: clean
clean:
	rm -f *png *dat *so

.PHONY: help
help:
	@echo "For comparison between e.g., Brian2 and Arbor, run: make comparison_brian2_arbor.png"
