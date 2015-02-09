TARGET = ../../darwin/oxs
SRCS = $(shell ls *.cc)
HEADS = $(shell ls *.h)
OBJSDIR = ../../darwin
OBJS = $(addprefix $(OBJSDIR)/,$(subst .cc,.o,$(SRCS)))

LOCALDIR = ../../local
OOMMFDIR = ../../../..

TESTMIF = testdata/Hyz_02_run.mif

.SUFFIXES: .cc .h .o

RM = rm
CP = cp
OOMMF = tclsh $(OOMMFDIR)/oommf.tcl

all: $(OBJS)
	$(OOMMF) pimake -cwd $(OOMMFDIR)

$(OBJSDIR)/%.o: %.cc %.h 
	$(CP) $*.cc $(LOCALDIR)
	$(CP) $*.h $(LOCALDIR)

clean:
	$(RM) -f $(OBJS) $(TARGET) *~ \#*\#

test:
	$(OOMMF) oxsii $(TESTMIF)

# Dependency rule
#-include Makefile.depend
