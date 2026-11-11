function plotComparison(results,outputDir)
%PLOTCOMPARISON Native plots without Python, LaTeX, or Symbolic Toolbox.
modes={'bCBF','blending','OI'}; colors=lines(3); sample=results.OI;
p=oicbf.config(sample.model,sample.profile);
fig=figure('Visible','off','Position',[100 100 1100 260*(2+p.m)]); tiledlayout(2+p.m,1);
nexttile; hold on;
for i=1:3, r=results.(modes{i}); plot(r.t,r.h,'Color',colors(i,:)); end
yline(0,'k--'); ylabel('h(x)'); legend(modes,'Location','best'); grid on;
nexttile; hold on;
for i=1:3, r=results.(modes{i}); plot(r.t(1:end-1),r.mu,'Color',colors(i,:)); end
ylabel('\mu'); grid on;
for j=1:p.m
    nexttile; hold on;
    for i=1:3
        r=results.(modes{i}); plot(r.t(1:end-1),r.u(j,2:end),'Color',colors(i,:));
    end
    yline(p.lb(j),'k:'); yline(p.ub(j),'k:');
    ylabel(sprintf('u_%d',j)); grid on;
end
xlabel('time (s)');
exportgraphics(fig,fullfile(outputDir,'controls_and_safety.png')); close(fig);
fig=figure('Visible','off'); hold on;
for i=1:3
    r=results.(modes{i});
    if strcmp(r.model,'double_integrator'), plot(r.x(1,:),r.x(2,:),'Color',colors(i,:));
    else, plot(r.x(5,:)/1000,r.x(4,:)/1000,'Color',colors(i,:)); end
end
if strcmp(r.model,'double_integrator'), xline(0,'k--'); xlabel('position'); ylabel('velocity');
else, yline(25,'k--'); xlabel('east (km)'); ylabel('north (km)'); end
legend(modes,'Location','best'); grid on;
exportgraphics(fig,fullfile(outputDir,'trajectories.png')); close(fig);
fig=figure('Visible','off','Position',[100 100 1100 230*ceil(p.n/2)]);
tiledlayout(ceil(p.n/2),2);
if p.n==2
    labels={'position','velocity'};
else
    labels={'\phi (rad)','\theta (rad)','\psi (rad)','north (m)', ...
        'east (m)','altitude (m)','P (rad/s)','N_z'};
end
for j=1:p.n
    nexttile; hold on;
    for i=1:3, r=results.(modes{i}); plot(r.t,r.x(j,:),'Color',colors(i,:)); end
    ylabel(labels{j}); xlabel('time (s)'); grid on;
end
exportgraphics(fig,fullfile(outputDir,'state_histories.png')); close(fig);
end
