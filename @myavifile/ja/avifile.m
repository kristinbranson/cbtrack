%AVIFILE  �V���� AVI �t�@�C���̍쐬
%
%   AVIOBJ = AVIFILE(FILENAME) �́A�f�t�H���g�̃p�����[�^�l������ AVIFILE 
%   �I�u�W�F�N�g AVIOBJ ���쐬���܂��BFILENAME ���g���q���܂�ł��Ȃ��ꍇ�A
%   '.avi' ���g���܂��BAVIFILE �ŊJ���ꂽ�t�@�C�������ɂ́AAVIFILE/CLOSE ��
%   �g�p���܂��B�J���Ă��邷�ׂĂ� AVI �t�@�C�������ɂ́A"clear mex" ���g�p
%   ���Ă��������B
%
%   AVIOBJ = AVIFILE(FILENAME,'PropertyName',VALUE,'PropertyName',VALUE,...) 
%   �́A�w�肵���v���p�e�B�l������ AVIFILE �I�u�W�F�N�g��Ԃ��܂��B
%
%   AVIFILE �p�����[�^
%
%   FPS         - AVI ���[�r�[�p�̕b���̃t���[���B���̃p�����[�^�́AADDFRAME 
%   �̎g�p�O�ɐݒ肳��Ă��Ȃ���΂Ȃ�܂���B�f�t�H���g�́A15 fps �ł��B
%
%   COMPRESSION - ���k�ɗp������@���w�肷�镶����ł��BUNIX �ł́A���̒l�� 
%   'None' �łȂ���΂Ȃ�܂���BWindows �p�̗��p�\�ȃp�����[�^�́A
%   'Indeo3', 'Indeo5', 'Cinepak', 'MSVC', 'RLE' �܂��� 'None' �ł��B
%
%   ���[�U�̈��k�@���g�p����ꍇ�A�l�͕������̃h�L�������g�Ŏw�肳�ꂽ 4 ��
%   �����R�[�h���g�p���܂��B�w�肳�ꂽ���[�U�̈��k�@��������Ȃ��ꍇ�A
%   ADDFRAME ���Ăяo���Ă���ԂɃG���[�ɂȂ�܂��B���̃p�����[�^�́AADDFRAME ��
%   �g�p�O�ɐݒ肳��Ă��Ȃ���΂Ȃ�܂���B
%   �f�t�H���g�́AWindows �ł� 'Indeo5'�AUNIX �ł� 'None' �ł��B
%   ����: Indeo5 �́AWindows �̂������̃o�[�W�����ŗ��p�ł��Ȃ��\��������܂��B
%
%   QUALITY      - 0 ���� 100 �̊Ԃ̐��B���̃p�����[�^�͈��k����Ă��Ȃ�
%   ���[�r�[�ɂ͉e�����܂���B���̃p�����[�^�́AADDFRAME �̎g�p�O�ɐݒ肳���
%   ���Ȃ���΂Ȃ�܂���B�����l�́A�����r�f�I�掿�ł��傫�ȃt�@�C���T�C�Y�ɂȂ�A
%   �Ⴂ�l�́A�Ⴂ�r�f�I�掿�ł�菬���ȃt�@�C���T�C�Y�ɂȂ�܂��B
%   �f�t�H���g�� 75 �ł��B
%
%   KEYFRAME     - �ꎞ�I�Ȉ��k���T�|�[�g���鈳�k�@�ɑ΂��āA���̃p�����[�^
%   �́A�P�ʎ��� (�b) ������̃L�[�t���[�����ł��B���̃p�����[�^�́AADDFRAME ��
%   �g�p�O�ɐݒ肳��Ă��Ȃ���΂Ȃ�܂���B�f�t�H���g�͖��b 2 �t���[���ł��B
%
%   COLORMAP     - �C���f�b�N�X�t�� AVI ���[�r�[�ɑ΂��Ďg�p�����J���[�}�b�v
%   ���`���� M �s 3 ��̍s��ł��BM �� 256 ��菬�������łȂ���΂Ȃ�܂��� 
%   (Indeo ���k�� 236 �ł�)�B���̃p�����[�^�́AMATLAB �̃��[�r�[�V���^�b�N�X��
%   ���� ADDFRAME ���g���ȊO�ł́AADDFRAME ���Ăяo�����O�ɐݒ肳��Ă��Ȃ����
%   �Ȃ�܂���B�f�t�H���g�̃J���[�}�b�v�͂���܂���B
%
%   VIDEONAME    - �r�f�I�X�g���[���p�̋L�q���B���̃p�����[�^�́A64 �������
%   �������AADDFRAME ���g�p����O�ɐݒ肳��Ă��Ȃ���΂Ȃ�܂���B�f�t�H���g��
%   �t�@�C�����ł��B
%
%
%   AVIFILE �v���p�e�B�́AMATLAB �\���̃V���^�b�N�X��p���Đݒ�ł��܂��B
%   ���Ƃ��΁A�ȉ��̃V���^�b�N�X���g���� Quality �v���p�e�B�� 100 ��ݒ肵�܂��B
%
%      aviobj = avifile(filename);
%      aviobj.Quality = 100;
%
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
%
%   �Q�l AVIFILE/ADDFRAME, AVIFILE/CLOSE, MOVIE2AVI.


%   Copyright 1984-2008 The MathWorks, Inc.
