# Possible build files:
# ulxlcmod.c <-- simulation code
# ulxlc_run.c <-- single use ulxlc runner.

# Use the following to compile shared c files that can be run in python :DD
# gcc -shared -Wl,-soname,ulxlc -o ulxlc_shared.so -fPIC ulxlcbase.c ulxlc_shared.c -L /home/x1/cfitsio/lib -lcfitsio -lm -lnsl


ulxlcmodxspec: ulxlc_run.c ulxlcbase.c
	gcc -g ulxlcbase.c ulxlc_run.c -L /home/x1/cfitsio/lib -lcfitsio -lm -lnsl -l sqlite3
