% script to generate a test image and store it in tiff format...

x = linspace (0,1,11);
y1 = x.^2;
y2 = 1-y1;
y3 = 4*x.*(1-x);
y4 = 1-y3;

plot (x,y1,'k:+');
hold on;
plot (x,y2,'k:x');
plot (x,y3,'k:o');
plot (x,y4,'k:*');

text (x(6)+0.03,y1(6),'{\it{y}}={\it{x}}^2');
text (x(6)+0.03,y2(6),'{\it{y}}=1-{\it{x}}^2');
text (x(3)+0.03,y3(3),'{\it{y}}=4{\it{x}}(1-{\it{x}})');
text (x(3)+0.03,y4(3),'{\it{y}}=1-4{\it{x}}(1-{\it{x}})');

xlabel ('{\it{x}}');
ylabel ('{\it{y}}');

print -dtiff test1
