NLP Programming Assignment #1 Write-up

Question 4:

The function rare_replace() performs the task of replacing words in Hmm's emission_counts member. It looks for a given word across multiple tags when considering if the word is rare. If across all tags the total is less than 5 occurrences it is deleted from the dictionary. These new counts are dumped to a file called 'new_counts'.

The function baseline_tag() runs through the keys of emiss_prob_log which stores the emission probabilities and grabs the highest scoring word/tag probability. If there is a word that isn't in the new dictionary i'm constructing (called emiss_high_prob) or if this word/tag probability is higher than the current one in emiss_high_prob it is added. Now we have
the highest probability for for each word and its corresponding tag.

The function emission_gen() goes through all of the counts in Hmm.emission_counts and logs the emission probabilities.

The function write_baseline_tagger() runs through the development data and in the event that a word is no in emiss_high_prob, the word and the "_rare_" probability are output. Otherwise the normal word and probability are output from emiss_high_prob. The output is a file called 'baseline_output'.

Baseline entity tagger evaluation:

     		
  precision     recall    	     F1-Score
Total:     	0.221961    0.525544    0.312106
PER:     	0.435451    0.231230    0.302061
ORG:     	0.475936    0.399103    0.434146
LOC:     	0.147750    0.870229    0.252612
MISC:     	0.491689    0.610206    0.544574

As a baseline this gives a good indication of how our tagger should improve once we 
introduce trigram context. As it was stated on Piazza we're using only emission probabilities to predict tags, so unfortunately we're choosing the same tag for every rare word.

Problem 5:

Using the updated counts  the function trigram_file_est() produces trigram estimates by ingesting 'ner_counts.dat' and ignoring all lines that don't contain "3-gram" as the second word. 

The function trigram_estimate() is a helper function that actually returns the probability estimate. The trigram estimates are dumped to a file called 'trigram_ests'.

These results of using the rare tagging is as follows:

Found 4704 NEs. Expected 5931 NEs; Correct: 3647.

	 precision 	recall 		F1-Score
Total:	 0.775298	0.614905	0.685849
PER:	 0.762535	0.595756	0.668907
ORG:	 0.611855	0.478326	0.536913
LOC:	 0.876458	0.696292	0.776056
MISC:	 0.830065	0.689468	0.753262

Overall the tagger's performed much better than the naive tagger. It is still lacking in performance in tagging ORG's.

This is the performance after introducing categories 'capInit' for names and 'num' for numbers. 

Found 5586 NEs. Expected 5931 NEs; Correct: 4230.

	 precision 	recall 		F1-Score
Total:	 0.757250	0.713202	0.734566
PER:	 0.806140	0.771491	0.788435
ORG:	 0.563624	0.632287	0.595985
LOC:	 0.844148	0.723555	0.779213
MISC:	 0.847480	0.693811	0.762985

Although total precision went down there was a significant improvement to recall for ORG's. Total F1 went up as well.

I left these two categories in after finding that other things such as categories for dates, All-caps, etc. weakened the performance significantly. Maybe they're not prevalent enough to make a major difference.

There are two shell scripts to run:

./problem4_5.sh runs numbers 4 and 5. Two result files are generated - baseline_output is the output from question 4 and vit_output is the results from part 5.

./problem6.sh runs number 6 and the results dump to vit_output.

performs the requirements for
