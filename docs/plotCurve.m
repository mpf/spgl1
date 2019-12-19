f = @(x) max(0, (x-2.8).^2 - 0.2);
g = @(x) min(0, 2*(x-2.8));

xvals = [0.0, 0.0, 0.0, 1.2, 1.2, 1.2];
delta = [0.0, 0.0, 0.0, 0.0, 0.5,-0.5];
mult  = [1.0, 1.0, 1.0, 1.0, 0.9, 1.1];

% Two steps of root finding
for i = [1,2]
   x = xvals(i);
   y = f(x);
   v = g(x);
   xvals(i+1) = (sigma - y) / v + x;
end

tau = linspace(0,3,1001);

% Root finding
for i = 1:length(xvals)
   h = plot(tau,f(tau),'b-');
   set(h,'linewidth',2);
   set(gca,'Fontsize',16);
   xlabel('\tau'); ylabel('Residual norm');
   ylim([0,10]);

   hold on;

   sigma = 1.0;
   h = plot([0,3],[sigma,sigma],'r--');
   set(h,'linewidth',2)

   x = xvals(i);
   y = f(x) + delta(i);
   v = g(x) * mult(i);
   
   tauBar = (sigma - y) / v + x;

   h = plot(x,y,'bo'); set(h,'MarkerSize',5,'linewidth',3,'color',0.7*[1,1,1]);
   h = plot([x,x],[0,10],'k--'); set(h,'color',0.7*[1,1,1]);
   h = plot(tau,y+(tau-x) * v,'k-');
   set(h,'linewidth',2,'color',0.7*[1,1,1]);
   h = plot(tauBar, sigma,'ro'); set(h,'MarkerSize',5,'linewidth',3);

   hold off;

   filename = sprintf('RootFind%d.png',i);
   exportFigure(filename);
end


% Duality gap
clf; gca; hold on; box on;
h = plot(tau,f(tau),'b-');
set(h,'linewidth',2);
set(gca,'Fontsize',16);
xlabel('\tau'); ylabel('Residual norm');
ylim([0,10]);

x = 1.5;
h = plot([x,x],[0,10],'k--'); set(h,'color',0.7*[1,1,1]);
h = plot(x, f(x)+0.2, 'm*', x, f(x)-0.4, 'm*');
set(h,'Markersize',6,'linewidth',2)

hold off
exportFigure('Lasso.png');


% Pareto curve
clf; gca; hold on; box on;
h = plot(tau,f(tau),'b-');
set(h,'linewidth',2);
set(gca,'Fontsize',16);
xlabel('\tau'); ylabel('Residual norm');
ylim([0,10]);
hold off
exportFigure('ParetoCurve.png');

hold on;
sigma = 1.0;
h = plot([0,3],[sigma,sigma],'r--'); set(h,'linewidth',2)

% Twenty steps of root finding
x = 0;
for i = 1:20
   y = f(x);
   v = g(x);
   x = (sigma - y) / v + x;
end
h = plot([x,x],[0,f(x)],'r-'); set(h,'Linewidth',2);
h = plot(x,f(x),'r*'); set(h,'Markersize',6,'linewidth',2);
hold off;
exportFigure('ParetoSigma.png');




function exportFigure(filename)
   print(gcf, '-dpng',filename);
   img = imread(filename);
   
   colIdx = find(min(sum(img,3),[],1) ~= 3*255);
   rowIdx = find(min(sum(img,3),[],2) ~= 3*255);
   img = img(rowIdx(1)-1:rowIdx(end)+1, colIdx(1)-1:colIdx(end)+1, :);
   imwrite(img,filename);
end