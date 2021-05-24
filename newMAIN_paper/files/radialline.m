function radialline(theta,r)

r=linspace(0,r,r*100);
[x,y]=pol2cart(theta,r);


plot(x,y,'k')
end