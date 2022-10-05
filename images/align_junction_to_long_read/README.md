## Description

In SNAF, I basically ask the question, given a local junction event, which long-read transcript it can correspond to.
Although in the script, I always work with AltAnalyze splicing ID (i.e. ENSG00001:E4.5-E6.2), the concept should be applied to
any junction represented by its chromsome coordiates (i.e. chr1:666666-777777).

The function is [recover_full_length_protein_long_read](https://github.com/frankligy/SNAF/blob/439f860ac123b68c89648624c7ac9788b9d91b84/snaf/surface/main.py#L1162-L1231). As you can see at the first glance, this function takes an argument called `gtf_dict`, which is a dictionary built from a standard GTF file
coming off from PacBio post-SQUANTi results. We will get to how this is built later, but the purpose is to store useful information (exon coordinates, gene id, etc) in a dictionary so that retrieval can be fast.

The logic of this function is basically, given the junction coordinate, we first narrow the potentially matched long-read transcript to the ones residing on 
same chromsome and strand with the junction, see [this line](https://github.com/frankligy/SNAF/blob/439f860ac123b68c89648624c7ac9788b9d91b84/snaf/surface/main.py#L1173). Then we basically iterate all the possible long-read transcript, and ask whether the junction coordinate match up with any ends of the exons in long-read gtf. As a step 1, I further narrowed down the potential transcript by eliminating the ones whose range is just out of short read junction, say the long-read is from 100-200, whereas the tested junction is from 400-500. When we have all `candidate_transcripts`, we start to see how they match, imagine one long-read transcript whose exons coordiantese are [(10,20),(30,40),(50,60)], if a short read junction is 40-50, then it will match to this long-read junction, so we basically make an array of the exons coordinate of each candidate transcript, get the index of start and end coordiante of short read junction, and see if they are adjacent or not (ei = si + 1).

The above logic doesn't solve all the scenarios, imaging if you have intron retention (ENSG00001:E4.5-I4.1, meaning the Intron4 is retained in the transcript), now when we consider the short read junction coordinates, we should actually consider the start as the start of E4.1 (the start of current exon), and the end is the end of E5.1 (the next exon instead of the intron itself). Because if the long-read transcript is [(10,20),(30,40),(50,60)] again, and the intron retention event is retaining the intron2, so the short read junction we consider should be 30-60, instead of 40-50. It is a bit confusing but this is how I do that in my script, not ideal but it works.

Now the remaining question, how to build the `gtf_dict`, it is handled by function [process_est_or_long_read_with_id](https://github.com/frankligy/SNAF/blob/439f860ac123b68c89648624c7ac9788b9d91b84/snaf/surface/main.py#L625).

I have some other related scripts, but I think the above is the most relevant one. 