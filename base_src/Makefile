# Environment-dependent directories
TARGET = ../../../darwin/oxs
OBJSDIR = ../../../darwin

BASEDIR = ../../../base
KEEPORIGSRCSDIR = ./original
OOMMFDIR = ../../../../..

SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)

OBJS = $(addprefix $(OBJSDIR)/,$(subst .cc,.o,$(SRCS)))

.SUFFIXES: .cc .h .o

RM = rm
CP = cp
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

.PHONY: all clean cleantest
all: $(OBJS)

$(OBJSDIR)/%.o: %.cc %.h 
	@if [ ! -d $(KEEPORIGSRCSDIR) ]; then mkdir $(KEEPORIGSRCSDIR); fi
	@if [ ! -f $(KEEPORIGSRCSDIR)/$*.cc ]; then \
		echo Copying the original file $*.cc to $(KEEPORIGSRCSDIR).; \
		$(CP) $(BASEDIR)/$*.cc $(KEEPORIGSRCSDIR)/; \
	fi
	@if [ ! -f $(KEEPORIGSRCSDIR)/$*.h ]; then \
		echo Copying the original file $*.h to $(KEEPORIGSRCSDIR).; \
		$(CP) $(BASEDIR)/$*.h $(KEEPORIGSRCSDIR)/; \
	fi
	$(foreach i,$?,$(CP) $(i) $(BASEDIR);)

clean:
	$(RM) -f $(OBJS) *~
	$(CP) -p $(KEEPORIGSRCSDIR)/* $(BASEDIR)/

cleantest:
	@echo cleantest base_src

