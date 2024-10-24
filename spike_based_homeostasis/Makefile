.DEFAULT_GOAL := init

CONFIG := config.json
NUM_TRIALS := 50

init:
ifeq ($(strip $(CASE)),target)
	@echo "Running target case..."
	@$(MAKE) homeostasis.svg
else ifeq ($(strip $(CASE)),control)
	@echo "Running control case..."
	@$(MAKE) no_homeostasis.svg
else
	@echo "Running target and control case..."
	@$(MAKE) homeostasis.svg no_homeostasis.svg
endif

homeostasis-catalogue.so: $(wildcard mechanisms/*.mod)
	./build_arbor_catalogue.sh homeostasis mechanisms

arbor_input_0.dat arbor_traces_0.dat arbor_spikes_0.dat: homeostasis-catalogue.so $(CONFIG) arbor_homeostasis.py
	./run_arbor.sh $(CONFIG) $(NUM_TRIALS) 0 --catalogue ./$<

arbor_input_1.dat arbor_traces_1.dat arbor_spikes_1.dat: homeostasis-catalogue.so $(CONFIG) arbor_homeostasis.py
	./run_arbor.sh $(CONFIG) $(NUM_TRIALS) 1 --catalogue ./$<

brian2_input_0.dat brian2_traces_0.dat brian2_spikes_0.dat: $(CONFIG) brian2_homeostasis.py
	./run_brian2.sh $(CONFIG) $(NUM_TRIALS) 0

brian2_input_1.dat brian2_traces_1.dat brian2_spikes_1.dat: $(CONFIG) brian2_homeostasis.py
	./run_brian2.sh $(CONFIG) $(NUM_TRIALS) 1

no_homeostasis.svg: arbor_input_0.dat arbor_traces_0.dat arbor_spikes_0.dat brian2_input_0.dat brian2_traces_0.dat brian2_spikes_0.dat compare.py
	./compare.py $(CONFIG) $(NUM_TRIALS) 0 --save no_homeostasis.svg

homeostasis.svg: arbor_input_1.dat arbor_traces_1.dat arbor_spikes_1.dat brian2_input_1.dat brian2_traces_1.dat brian2_spikes_1.dat compare.py
	./compare.py $(CONFIG) $(NUM_TRIALS) 1 --save homeostasis.svg

.PHONY: init clean
clean:
	rm -f *svg *so *dat
