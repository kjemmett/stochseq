function plot_sequence_inference_static(o,m)
% function inf_obj = plot_sequence_inference_static(est,dna)
% 
% Display image object of inferred sequence against
% true DNA sequence
%
% Can be used independently or as a subroutine 
% for animation generation
% 
%
% inputs:
%   est : sequence estimate (Lx4)
%   dna : true dna sequence

est = o.S{end};
dna = m.dna;

figure;


% true sequence
ax(1)=subplot(3,1,1);
dna_mat = gen_dnamatrix(dna);
dna_obj = imagesc(dna_mat');
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Nucleotide');
set(gca, 'YTick', [1 2 3 4], 'YTickLabel',{'A','C','G','T'});

% inferred sequence
ax(2)=subplot(3,1,2);
inf_obj=imagesc(est');
hold on;
plot(dna','y-','LineWidth',2)
hold off;
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Nucleotide');
set(gca,'YTick',[1 2 3 4],'YTickLabel',{'A','C','G','T'});

% inference entropy
ax(3)=subplot(3,1,3);
ent_obj=bar(calc_entropy(est));
xlim([0 length(dna)]);
ylim([0 1]);
set(gca, 'FontSize', 14);
xlabel('Sequence Position');
ylabel('Entropy');

set(zoom,'Motion','horizontal','Enable','on')
linkaxes([ax(1) ax(2) ax(3)],'x');

set(gcf, 'Position', [100 100 1000 400]);
