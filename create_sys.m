 
load('sys_IO','f11');
[b a] = ss2tf(f11.a,f11.b,f11.c,f11.d);
sys = tf(b,a);

save('sys_IO','f11');