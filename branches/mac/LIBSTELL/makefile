include ../makefile.include

SRCS = $(foreach dir, $(shell ls), $(wildcard $(dir)/*.f*))
OBJS = \
	$(foreach dir, $(shell ls Sources), $(filter-out %.f90, $(patsubst  Sources/$(dir)/%.f, %.o, $(wildcard  Sources/$(dir)/*.f*)))) \
	$(foreach dir, $(shell ls Sources), $(filter-out %.f, $(patsubst  Sources/$(dir)/%.f90, %.o, $(wildcard  Sources/$(dir)/*.f*))))
	
CURRENT_DIR_NAME = $(patsubst $(ROOT_DIR)/%,%, $(shell pwd))
	
all: $(BUILD_DIR)/$(OBJS)

include $(BUILD_DIR)/$(CURRENT_DIR_NAME).dep

$(OBJS): $(patsubst ./%, %, $(shell find Sources -name "$(patsubst %.o,%, $@).*"))
	@echo Building $@
	@$(FC) $(FFLAGS) -c $(patsubst ./%, %, $(shell find Sources -name "$(patsubst %.o,%, $@).*")) -o $(BUILD_DIR)/$@
	
clean:
	@echo Deleting $(CURRENT_DIR_NAME) build files
	@-rm $(BUILD_DIR)/*.o *.mod
