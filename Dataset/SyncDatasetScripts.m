
%From BigPurple to NYUShare: direct

rsync -auvzPLK --exclude=WaveSpec --exclude=Kilosort* --exclude=GitHub /mnt/BigPurpleWilliam/ /mnt/NyuShare/dl2820/WMDataset/

rsync -auvzPLK --exclude=WaveSpec --exclude=*.dat --exclude=Kilosort*  /home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/OnCluster/ /home/dlevenstein/ProjectRepos/ACh-and-CorticalState/Dataset/OnDesktop/
