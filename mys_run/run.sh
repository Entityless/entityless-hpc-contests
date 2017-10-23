set -x
rm run.log
cp ../mys_code/yooooooo ./
bsub -I -J mys -o run.log -q q_sw_expr -n 1 -np 1 -share_size 4096 -host_stack 2048 -b -cgsp 64 ./yooooooo
