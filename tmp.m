x=[2;3]; b1=[1;0]; b2=[0;1];
figure();
subplot(121);
hold on;č

% promitnuti do bazi a vypocteni promitnutych souradnic
x1=x'*b1
x2=x'*b2
cx1=b1.*x1;
cx2=b2.*x2;

p1=plot(cx1(1),cx1(2),'MarkerFaceColor',[0 0 1],'Marker','square', 'Color',[0 0 1]);
set(p1,'Clipping','off')
p2=plot(cx2(1),cx2(2),'MarkerFaceColor',[0 0 1],'Marker','square', 'Color',[0 0 1]);
set(p2,'Clipping','off')
daspect([1 1 1]); % zachovani pomeru vsech os v plotu

% do dalsiho grafu vektor z promitnutych souradnic zesyntetizujeme
subplot(122);
hold on;

% secteme promitnute souradnice a tak dostaneme puvodni vektor x
synt_x=cx1+cx2;
plotv(synt_x,'g-'); % zelena
daspect([1 1 1]); % zachovani pomeru vsech os v plotu
hold off;