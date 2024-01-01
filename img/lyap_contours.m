figure(2)
clf

t = tiledlayout(1, 3)
nexttile
hold on
for i = 1:10
    plot([i, 0], [0, i], 'k')    
end
title('$x_1+x_2$', 'interpreter', 'latex')
pbaspect([1,1,1])
axis off

nexttile 
hold on
for i = 1:10
    plot([i, i, 0], [0, i, i], 'k')    
end
title('$\max(x_1,x_2)$', 'interpreter', 'latex')
axis off
pbaspect([1,1,1])

nexttile
hold on
th = linspace(0, pi/2, 100);
cx = cos(th);
cy = sin(th);
for i = 1:10
    plot(cx*i, cy*i, 'k')
end
title('$\sqrt{x_1^2 + x_2^2}$', 'interpreter', 'latex')
axis off
pbaspect([1,1,1])