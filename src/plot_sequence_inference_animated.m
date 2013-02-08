function plot_sequence_inference_animated(inf,model)
% function plot_sequence_inference_animated(inf,model)
% 
% Animate sequence inference
%
% inputs:
%   inf : struct containing inference output
%       inf.h : struct containing S estimate and entropy at each iteration
%   model : struct containing model parameters
%       model.dna : input DNA sequence

L = model.seqlength;
dna = model.dna;

figure
subplot(2,1,1);inf_plot=imagesc(0.25*ones(4,L));hold on;plot(dna','y.-');hold off;xlabel('Sequence Position');ylabel('Nucleotide');set(gca,'YTick',[1 2 3 4],'YTickLabel',{'A','C','G','T'});
subplot(2,1,2);ent_plot=bar(ones(1,L));xlim([0 length(dna)]);ylim([0 1]);xlabel('Sequence Position');ylabel('Entropy');

set(zoom,'Motion','horizontal','Enable','on')

pause

for i=1:size(inf.h,2)
	set(inf_plot,'CData',inf.h(i).S');
    %tot_ent = sum(inf.h(i).inf_ent);
	%title(['iteration ',num2str(i),' total entropy ',num2str(tot_ent)]);
	%set(ent_plot,'ydata',inf.h(i).inf_ent');
    pause
end
