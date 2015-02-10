TARGET = ../../darwin/oxs
SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
OBJSDIR = ../../darwin
SUBDIRS = ./base_src ./ext_src
LOCALDIR = ../../local
OOMMFDIR = ../../../..

SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
BASESRCS = $(shell ls $(BASESRCSDIR)/*.cc)

OBJS = $(addprefix $(OBJSDIR)/,$(subst .cc,.o,$(SRCS)))

TESTMIF = test01.mif

.SUFFIXES: .cc .h .o

RM = rm
CP = cp
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

.PHONY: all clean test $(SUBDIRS)
all: $(SUBDIRS) $(OBJS)
	$(OOMMF) pimake -cwd $(OOMMFDIR)

$(SUBDIRS):
	@echo Running make in subdirectory $@
	@$(MAKE) -C $@

$(OBJSDIR)/%.o: %.cc %.h 
	$(CP) $*.cc $(LOCALDIR)
	$(CP) $*.h $(LOCALDIR)

clean:
	$(foreach i,$(SUBDIRS),$(MAKE) -C $(i) clean;)
	$(RM) -f $(OBJS) $(TARGET) *~
	$(RM) -f $(addprefix $(LOCALDIR)/,$(SRCS))
	$(RM) -f $(addprefix $(LOCALDIR)/,$(HEADS))

test:
	$(OOMMF) oxsii $(TESTMIF)

