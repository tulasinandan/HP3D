PARAMFILE=param_hybrid
pex=`awk '/^#define/ && / pex / {print $3}' $PARAMFILE`
pey=`awk '/^#define/ && / pey / {print $3}' $PARAMFILE`
nx=`awk '/^#define/ && / nx / {print $3}' $PARAMFILE`
ny=`awk '/^#define/ && / ny / {print $3}' $PARAMFILE`
lx=`awk '/^#define/ && / lx / {print $3}' $PARAMFILE`
ly=`awk '/^#define/ && / ly / {print $3}' $PARAMFILE`
movieout=`awk '/^#define/ && / movieout / {print $3}' $PARAMFILE`
movieout_full=`awk '/^#define/ && / movieout_full / {print $3}' $PARAMFILE`
dt=`awk '/^#define/ && / dt / {print $3}' $PARAMFILE`
d_e2=`awk '/^#define/ && / d_e2 / {print $3}' $PARAMFILE`
d_b=`awk '/^#define/ && / D_b / {print $3}' $PARAMFILE`
d_pe=`awk '/^#define/ && / D_pe / {print $3}' $PARAMFILE`
T_e=`awk '/^#define/ && / T_e / {print $3}' $PARAMFILE`
T_i=`awk '/^#define/ && / T_i / {print $3}' $PARAMFILE`
ppg=`awk '/^#define/ && / ppg / {print $3}' $PARAMFILE`

