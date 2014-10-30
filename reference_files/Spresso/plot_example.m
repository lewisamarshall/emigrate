load Input' Files'/Case1_ITPDemo.mat
t=[15, 50, 100];
mu=[30,10,0,20];
for i=1:3
    subplot(3,3,i)
    hold on
    for j=[1,2,3]
        plot(cMatAllTimes(:,j,t(i)));
    end
    ylim([0, 15])
end

for i=1:3
    subplot(3,3,3+i)
    plot(SigVecAllTime(:,t(i)));
end

for i=1:3
    subplot(3,3,6+i)
    hold on
    for j=[1,2,4]
        plot(mu(j)./SigVecAllTime(:,t(i)));
    end
end
    