from dna_features_viewer import BiopythonTranslator
import os 

for file in os.listdir('dashboard/static/annotations'):
	print(file.split('.')[0])
	graphic_record = BiopythonTranslator().translate_record('dashboard/static/annotations/'+file)
	ax, _ = graphic_record.plot(figure_width=30, strand_in_label_threshold=7)

	ax.figure.savefig('dashboard/static/annotations/'+file.split('.')[0]+'.png', bbox_inches='tight')

