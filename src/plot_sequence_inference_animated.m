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

fighandle = figure('Position', [100 100 1000 400]);
set(fighandle, 'PaperPositionMode', 'auto')

% true sequence
ax(1)=subplot(3,1,1);
dna_mat = gen_dnamatrix(model.dna);
dna_obj = imagesc(dna_mat');
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Nucleotide');
set(gca, 'YTick', [1 2 3 4], 'YTickLabel', {'A','C','G','T'});

% inferred sequence
ax(2)=subplot(3,1,2);
inf_plot=imagesc(0.25*ones(4,L));
%hold on;
%plot(dna','y.-');
%hold off;
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Nucleotide');
set(gca,'YTick',[1 2 3 4],'YTickLabel',{'A','C','G','T'});

% inference entropy
ax(3)=subplot(3,1,3);
ent_plot=bar(ones(1,L));
xlim([0 length(dna)]);
ylim([0 1]);
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Entropy');

set(zoom,'Motion','horizontal','Enable','on')
linkaxes([ax(1) ax(2) ax(3)], 'x');

print -dpng '/Users/kje/Desktop/test/test00.png'

for i=1:size(inf.S,2)
	set(inf_plot,'CData',inf.S{i}');
	%title(['iteration ',num2str(i),' total entropy ',num2str(tot_ent)]);
	set(ent_plot,'ydata',calc_entropy(inf.S{i}));
    %pause(.5)
    ofile = ['/Users/kje/Desktop/test/test' sprintf('%02d', i) '.png'];
    print('-dpng', ofile);
end

