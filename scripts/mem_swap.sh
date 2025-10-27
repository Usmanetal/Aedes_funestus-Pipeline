# sudo fallocate -l 8G /swapfile2

# # set permission
# sudo chmod 600 /swapfile2

# # make a swap space
# sudo mkswap /swapfile2

# # Enable it
# sudo swapon /swapfile2

# # keep watching:
# # watch -n 1 free -m

sudo fallocate -l 16G /swapfile2
sudo chmod 600 /swapfile2
sudo mkswap /swapfile2
sudo swapon /swapfile2

