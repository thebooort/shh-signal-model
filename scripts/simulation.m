x = lsode ("f", [0; 0;0;0.01], (t = linspace (0, 2000, 20000)'));
plot(t,x)