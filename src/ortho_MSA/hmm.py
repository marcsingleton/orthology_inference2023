"""Classes and functions for discriminative HMM training."""

from functools import reduce
from itertools import accumulate

import homomorph
from numpy import exp, log


class HMM(homomorph.HMM):
    """A class that modifies the default HMM to accept pre-calculated intermediates.

    Discriminative training requires decoding of both states and transitions.
    Both of these calculations rely on the forward and backward variables.
    Since these are the most computationally intensive steps in decoding and
    both are used in state and transition decoding, this modified class accepts
    the forward and backward variables as arguments rather than calculating
    them from scratch each time.
    """
    def forward_backward1(self, emits, fs, ss_f, bs, ss_b):
        """Posterior decoding of states."""
        p = reduce(lambda x, y: x+y, map(log, ss_f))
        ss_f = list(accumulate(map(log, ss_f)))
        ss_b = list(accumulate(map(log, ss_b[::-1])))[::-1]

        fbs = {state: [] for state in self.states}
        for i in range(len(emits)):
            for state, fb in fbs.items():
                fbs[state].append(fs[state][i]*bs[state][i]*exp(ss_f[i]+ss_b[i]-p))

        return fbs

    def forward_backward2(self, emits, fs, ss_f, bs, ss_b):
        """Posterior decoding of transitions."""
        p = reduce(lambda x, y: x+y, map(log, ss_f))
        ss_f = list(accumulate(map(log, ss_f)))
        ss_b = list(accumulate(map(log, ss_b[::-1])))[::-1]

        pairs = [(state0, state1) for state0, t_dist in self.t_dists.items() for state1 in t_dist]
        fbs = {pair: [] for pair in pairs}
        for i in range(len(emits)-1):
            for (state0, state1), fb in fbs.items():
                fb.append(fs[state0][i]*self.t_dists[state0][state1]*self.e_dists[state1].pmf(emits[i+1])*bs[state1][i+1]*exp(ss_f[i]+ss_b[i+1]-p))

        return {pair: sum(fb) for pair, fb in fbs.items()}

    def joint_likelihood(self, emits, states):
        """Log-likelihood of emission and state sequence."""
        states = [self._state2idx[state] for state in states]
        p, state0 = log(self._e_dists_rv[states[0]].pmf(emits[0]) * self._start_dist_rv.pmf(states[0])), states[0]
        for emit, state1 in zip(emits[1:], states[1:]):
            p += log(self._e_dists_rv[state1].pmf(emit) * self._t_dists_rv[state0].pmf(state1))
            state0 = state1
        return p


def count_transitions(state_seq, t_sets):
    """Return counts of transitions between states."""
    mijs = {(state0, state1): 0 for state0, t_set in t_sets.items() for state1 in t_set}
    state0 = state_seq[0]
    for state1 in state_seq[1:]:
        mijs[(state0, state1)] += 1
        state0 = state1
    return mijs


def count_states(state_seq, state_set):
    """Return counts of states."""
    mis = {state: [] for state in state_set}
    for state in state_seq:
        for s in state_set:
            mis[s].append(1 if state == s else 0)
    return mis
