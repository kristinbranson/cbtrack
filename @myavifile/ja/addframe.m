%ADDFRAME  �r�f�I�t���[���� AVI �t�@�C���ɒǉ�
%
%   AVIOBJ = ADDFRAME(AVIOBJ,FRAME) �́AFRAME ���̃f�[�^�� AVIFILE ��
%   �쐬���ꂽ AVIOBJ �ɒǉ����܂��BFRAME �̓C���f�b�N�X�t���C���[�W (M�~N)�A
%   �܂��́Adouble �� uint8 �̐��x�̃g�D���[�J���[�C���[�W (M�~N�~3) ��
%   �����ꂩ�ɂȂ�܂��BFRAME �� AVI �t�@�C���� 1 �Ԗڂ̃t���[���ɒǉ�����Ȃ�
%   �ꍇ�́A�O�̃t���[���̎����Ƃ̐�������ۂ��Ȃ���΂Ȃ�܂���B
%
%   AVIOBJ = ADDFRAME(AVIOBJ,FRAME1,FRAME2,FRAME3,...) �́A�����̃t���[���� 
%   AVI �t�@�C���ɒǉ����܂��B
%
%   AVIOBJ = ADDFRAME(AVIOBJ,MOV) �́AMATLAB ���[�r�[ MOV �Ɋ܂܂��t���[��
%   �� AVI �t�@�C���ɒǉ����܂��B�C���f�b�N�X�t���C���[�W�Ƃ��ăt���[����
%   �i�[���� MATLAB ���[�r�[�́A�J���[�}�b�v�����炩���ߐݒ肳��Ă��Ȃ��ꍇ�A
%   AVI �t�@�C���p�̃J���[�}�b�v�Ƃ��� 1 �Ԗڂ̃t���[���̃J���[�}�b�v���g�p���܂��B
%
%   AVIOBJ = ADDFRAME(AVIOBJ,H) �́AFigure�A�܂��́A���̃n���h�� H ����
%   �t���[�����L���v�`�����A���̃t���[���� AVI �t�@�C���ɒǉ����܂��B
%   �t���[���́AAVI �t�@�C���ɒǉ������O�ɁA��ʊO�̔z��ɕ`�悳��܂��B
%   �A�j���[�V�������̃O���t�B�b�N�X�� XOR �O���t�B�b�N�X���g�p���Ă���ꍇ�A
%   ���̃V���^�b�N�X���g�p����K�v�͂���܂���B
%
%   �A�j���[�V������ XOR �O���t�B�b�N�X���g�p���Ă���ꍇ�AMATLAB ���[�r�[�� 
%   1 �̃t���[���ɃO���t�B�b�N�X���L���v�`���������� GETFRAME ���g�p���A
%   ���L�̗�̂悤�ɁA[AVIOBJ] = ADDFRAME(AVIOBJ,MOV) ���g�p���܂��B
%   GETFRAME �́A��ʏ�̃C���[�W�̃X�i�b�v�V���b�g���B��܂��B
%
%   ��:
%
%      t = linspace(0,2.5*pi,40);
%      fact = 10*sin(t);
%      fig=figure;
%      aviobj = avifile('example.avi')
%      [x,y,z] = peaks;
%      for k=1:length(fact)
%          h = surf(x,y,fact(k)*z);
%          axis([-3 3 -3 3 -80 80])
%          axis off
%          caxis([-90 90])
%          F = getframe(fig);
%          aviobj = addframe(aviobj,F);
%      end
%      close(fig)
%      aviobj = close(aviobj);
%
%   �Q�l AVIFILE, AVIFILE/CLOSE, MOVIE2AVI.


%   Copyright 1984-2008 The MathWorks, Inc.
