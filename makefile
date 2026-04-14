objects=$(OBJDIR)/tbfor__var.o \
        $(OBJDIR)/tbfor_allo.o \
        $(OBJDIR)/tbfor_prnt.o \
	$(OBJDIR)/tbfor_lalg.o \
        $(OBJDIR)/tbfor_read.o \
        $(OBJDIR)/tbfor_kvec.o \
        $(OBJDIR)/tbfor_init.o \
        $(OBJDIR)/tbfor_rcov.o \
        $(OBJDIR)/tbfor_hmlt.o \
        $(OBJDIR)/tbfor_cdos.o \
        $(OBJDIR)/tbfor_task.o \
        $(OBJDIR)/tbfor_self.o \
        $(OBJDIR)/tbfor_hmat.o \
        $(OBJDIR)/tbfor_post.o \
        $(OBJDIR)/tbfor_join.o \
        $(OBJDIR)/tbfor___hs.o \
        $(OBJDIR)/tbfor___es.o \
        $(OBJDIR)/tbfor___ex.o \
        $(OBJDIR)/tbfor.o
SRCDIR=Src
OBJDIR=Obj
LAPACK= ./liblapack.a      
BLAS= ./librefblas.a   
FC=mpif90
FFLAGS= -g -pg -fbounds-check -fbackslash 
tbfor: $(objects)  
	$(FC) ${FFLAGS} -o  tbfor $(objects) $(LAPACK) $(BLAS) 
$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) ${FFLAGS} -J $(OBJDIR) -c -o $@ $< 
clean:
	rm tbfor $(objects) $(OBJDIR)/*.mod     
 
