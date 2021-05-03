gcc -shared -Wl,-soname,ulxlc -o ulxlc_shared.so -fPIC ulxlcbase.c ulxlc_shared.c -L /home/x1/cfitsio/lib -lcfitsio -lm -lnsl
