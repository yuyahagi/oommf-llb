TARGET = ../../darwin/oxs
SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
OBJSDIR = ../../darwin
SUBDIRS = ./base_src ./ext_src
LOCALDIR = ../../local
OOMMFDIR = ../../../..

SRCS = $(shell ls *.cc | grep -v '\-threaded.cc')
HEADS = $(shell ls *.h)
BASESRCS = $(shell ls $(BASESRCSDIR)/*.cc)

OBJS = $(addprefix $(OBJSDIR)/,$(subst .cc,.o,$(SRCS)))

TESTMIF01 = test01.mif
TESTMIF02 = test02.mif
TESTMIF03 = test03.mif

.SUFFIXES: .cc .h .o

RM = rm
CP = cp
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

.PHONY: all clean test $(SUBDIRS)
all: $(SUBDIRS) $(OBJS) $(OBJSDIR)/yy_2latdemag-threaded.o
	$(OOMMF) pimake -cwd $(OOMMFDIR)

$(SUBDIRS):
	@echo Running make in subdirectory $@
	@$(MAKE) -C $@

$(OBJSDIR)/%.o: %.cc %.h 
	$(CP) $*.cc $(LOCALDIR)
	$(CP) $*.h $(LOCALDIR)

$(OBJSDIR)/yy_2latdemag-threaded.o: yy_2latdemag-threaded.cc yy_2latdemag.h
	$(CP) yy_2latdemag-threaded.cc $(LOCALDIR)

clean:
	$(foreach i,$(SUBDIRS),$(MAKE) -C $(i) clean;)
	$(RM) -f $(OBJS) $(TARGET) *~
	$(RM) -f $(addprefix $(LOCALDIR)/,$(SRCS))
	$(RM) -f $(addprefix $(LOCALDIR)/,$(HEADS))

test:
	$(OOMMF) oxsii $(TESTMIF01)

test02:
	@echo 'test02'
	$(OOMMF) oxsii $(TESTMIF02)

test03:
	$(OOMMF) oxsii $(TESTMIF03)

-include $(DEPEND)
