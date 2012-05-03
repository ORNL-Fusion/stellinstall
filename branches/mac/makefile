include makefile.include

ROOT_DIR = $(shell pwd)
export ROOT_DIR

all: v3fit

v3fit: vmec
	@echo -------------------- Creating $@ --------------------
	@$(MAKE) -C V3FITA
	@$(FC) $(LDFLAGS) -lstell -lvmec -Lbin -o bin/$@ V3FITA/$(BUILD_DIR)/*.o
	@echo -------------------- Created $@ --------------------

vmec: libstell.a libvmec.a
	@echo -------------------- Creating $@ --------------------
	@$(FC) $(LDFLAGS) -lstell -lvmec -Lbin -o bin/$@ VMEC2000/$(BUILD_DIR)/vmec.o
	@mv VMEC2000/*.mod V3FITA
	@echo -------------------- Created $@ --------------------

libvmec.a: bin
	@echo -------------------- Creating $@ --------------------
	@$(MAKE) -C VMEC2000
	@ar rcs bin/$@ $(patsubst %.o, VMEC2000/$(BUILD_DIR)/%.o, $(filter-out vmec.o, $(foreach dir, $(shell ls VMEC2000/Sources), $(filter-out %.f90, $(patsubst  VMEC2000/Sources/$(dir)/%.f, %.o, $(wildcard  VMEC2000/Sources/$(dir)/*.f*)))) $(foreach dir, $(shell ls VMEC2000/Sources), $(filter-out %.f, $(patsubst  VMEC2000/Sources/$(dir)/%.f90, %.o, $(wildcard  VMEC2000/Sources/$(dir)/*.f*))))))
	@echo -------------------- Created $@ --------------------

libstell.a: bin
	@echo -------------------- Creating $@ --------------------
	@$(MAKE) -C LIBSTELL
	@ar rcs bin/$@ LIBSTELL/$(BUILD_DIR)/*.o
	@mv LIBSTELL/*.mod VMEC2000
	@echo -------------------- Created $@ --------------------

bin:
	@mkdir bin
	
clean:
	@-rm -rf bin
	@$(MAKE) clean -C V3FITA
	@$(MAKE) clean -C VMEC2000
	@$(MAKE) clean -C LIBSTELL
	