# growTurNet
A variation of Scholes et al's method to determine the Turing space of various network topologies. This has been modified to study growth.

To use this code, first use "New_DefineSystem" and update the network details/parameters to what you want to study/loop over. Once this is saved and run, perform the actual simulations with the "New_RunNetwork" script. You will want to update the savePath variable to wherever you want to send the results, however. Depending on how many parameters you want to loop over, the simulations could take some time. Please contact me directly if you have any questions (ckonow@brandeis.edu).

The code and most functions are based on (and in some cases directly originate with) the paper cited as:
Scholes, N. S., Schnoerr, D., Isalan, M., & Stumpf, M. P. H. (2019). A Comprehensive Network Atlas Reveals That Turing Patterns Are Common but Not Robust. Cell Systems, 9(3), 243-257.e4. https://doi.org/10.1016/j.cels.2019.07.007

The modifications I have made to incorporate growth are my own, but are based on the mathematical analysis of Madzvamuse et al, originally in this paper:
Madzvamuse, A., Gaffney, E. A., & Maini, P. K. (2010). Stability analysis of non-autonomous reaction-diffusion systems: the effects of growing domains. Journal of Mathematical Biology, 61(1), 133–164. https://doi.org/10.1007/s00285-009-0293-4

As good as Madzvamuse's work is, I have found the introduction to Klika et al's paper (which contains a summary of Madzvamuse et al) easier to read and understand (since it deals less with proving the properties). While I don't worry so much about the history dependence, the way he states the problem and explains the math is more straightforward, so this paper may be worth taking a look at too:
Klika, V., & Gaffney, E. A. (2017). History dependence and the continuum approximation breakdown: the impact of domain growth on Turing’s instability. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Science, 473(2199), 20160744. https://doi.org/10.1098/rspa.2016.0744
